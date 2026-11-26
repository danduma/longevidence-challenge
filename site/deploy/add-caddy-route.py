"""Run on the existing host after uploading dist/challenge to /opt/longevidence-challenge/challenge.
Narrowly adds the two Longevidence hosts' /challenge path; preserves all fallback redirects.
"""
from pathlib import Path
import datetime,subprocess
config=Path('/etc/caddy/Caddyfile')
text=config.read_text()
old='''longevidence.org:443, www.longevidence.org:443, app.longevidence.org:443, the.longevidence.app:443, www.longevipedia.net:443 {
\ttls {
\t\tprotocols tls1.2 tls1.3
\t}

\tredir https://longevipedia.net{uri} permanent
}'''
new='''longevidence.org:443, www.longevidence.org:443, app.longevidence.org:443, the.longevidence.app:443, www.longevipedia.net:443 {
\ttls {
\t\tprotocols tls1.2 tls1.3
\t}

\t@longevidenceChallenge {
\t\thost longevidence.org www.longevidence.org
\t\tpath /challenge /challenge/*
\t}
\thandle @longevidenceChallenge {
\t\troot * /opt/longevidence-challenge
\t\theader X-Content-Type-Options nosniff
\t\theader Referrer-Policy same-origin
\t\theader Content-Security-Policy "default-src 'self'; style-src 'self'; script-src 'self'; img-src 'self' data:; connect-src 'none'; base-uri 'none'; frame-ancestors 'none'"
\t\tfile_server
\t}
\thandle {
\t\tredir https://longevipedia.net{uri} permanent
\t}
}'''
if new in text:
 print('Challenge route already present; no config change.');raise SystemExit(0)
if text.count(old)!=1:raise SystemExit('Expected redirect block changed; inspect config before proceeding.')
if not Path('/opt/longevidence-challenge/challenge/index.html').is_file():raise SystemExit('Static files missing')
backup=config.with_name('Caddyfile.backup-challenge-'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
backup.write_text(text);backup.chmod(0o600)
config.write_text(text.replace(old,new,1))
check=subprocess.run(['caddy','validate','--config',str(config)],capture_output=True,text=True)
if check.returncode:
 config.write_text(text);print(check.stderr);raise SystemExit('Validation failed; original config restored.')
action='restart' if '\tadmin off' in text else 'reload'
reload=subprocess.run(['systemctl',action,'caddy'],capture_output=True,text=True)
if reload.returncode:
 config.write_text(text);subprocess.run(['systemctl',action,'caddy'],check=True);raise SystemExit('Activation failed; original config restored and reactivated.')
print('Caddy validated and activated using',action,'Backup:',backup)
