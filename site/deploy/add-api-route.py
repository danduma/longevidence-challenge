"""Run on the host. Adds the submission-service proxy inside the existing /challenge handle.
Backs up the Caddyfile, validates, and restarts (admin API is off); restores on any failure.
"""
from pathlib import Path
import datetime, subprocess
config = Path('/etc/caddy/Caddyfile')
text = config.read_text()
old = '''\t\theader Content-Security-Policy "default-src 'self'; style-src 'self'; script-src 'self'; img-src 'self' data:; connect-src 'none'; base-uri 'none'; frame-ancestors 'none'"
\t\tfile_server
\t}'''
new = '''\t\theader Content-Security-Policy "default-src 'self'; style-src 'self'; script-src 'self'; img-src 'self' data:; connect-src 'none'; base-uri 'none'; frame-ancestors 'none'"
\t\t@challengeDynamic path /challenge/api/* /challenge/c/* /challenge/contributions /challenge/contributions/ /challenge/submit /challenge/submit/
\t\trequest_body @challengeDynamic {
\t\t\tmax_size 128KB
\t\t}
\t\treverse_proxy @challengeDynamic 127.0.0.1:8790
\t\tfile_server
\t}'''
dynamic = '/challenge/review /challenge/review/ /challenge/activity /challenge/activity/ /challenge/credits /challenge/credits/'
legacy_new = new
new = new.replace('/challenge/submit/\n', '/challenge/submit/ '+dynamic+'\n')
if new in text:
    print('API route already present; no config change.'); raise SystemExit(0)
if legacy_new in text:
    old = legacy_new

if text.count(old) != 1:
    raise SystemExit('Expected challenge handle changed; inspect config before proceeding.')
backup = config.with_name('Caddyfile.backup-challenge-api-' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
backup.write_text(text); backup.chmod(0o600)
config.write_text(text.replace(old, new, 1))
check = subprocess.run(['caddy', 'validate', '--config', str(config), '--adapter', 'caddyfile'], capture_output=True, text=True)
if check.returncode:
    config.write_text(text); print(check.stderr[-2000:]); raise SystemExit('Validation failed; original config restored.')
action = 'restart' if '\tadmin off' in text else 'reload'
run = subprocess.run(['systemctl', action, 'caddy'], capture_output=True, text=True)
if run.returncode:
    config.write_text(text); subprocess.run(['systemctl', action, 'caddy'], check=True)
    raise SystemExit('Activation failed; original config restored and reactivated.')
print('Caddy validated and activated using', action, 'Backup:', backup)
