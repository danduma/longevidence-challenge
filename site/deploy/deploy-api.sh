#!/bin/sh
# Install or update the submission service on the existing Longevidence host, then add its Caddy route (idempotent).
# Usage: SSH_KEY=/path/to/key ./deploy/deploy-api.sh   (then ./deploy/deploy-static.sh for the pages)
set -eu
cd "$(dirname "$0")/.."
HOST=${HOST:?Set HOST to the deployment SSH destination}
SSH_KEY=${SSH_KEY:?Set SSH_KEY to your deployment key path}
SSH="ssh -o BatchMode=yes -o ConnectTimeout=10 -i $SSH_KEY"

python3 -m py_compile api/server.py api/review_store.py api/review_views.py
$SSH "$HOST" python3 - <<'PYBACKUP'
import datetime, pathlib, shutil, sqlite3
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
api=pathlib.Path('/opt/longevidence-challenge/api')
if api.exists(): shutil.copytree(api,api.with_name('api.backup-'+stamp))
database=pathlib.Path('/var/lib/longevidence-challenge/challenge.db')
if database.exists():
    target=database.with_name('challenge.backup-'+stamp+'.db')
    with sqlite3.connect(database) as source, sqlite3.connect(target) as destination:
        source.backup(destination)
    target.chmod(0o600)
print('API and consistent SQLite backups:',stamp)
PYBACKUP
$SSH "$HOST" 'id lcapi >/dev/null 2>&1 || useradd --system --no-create-home --shell /usr/sbin/nologin lcapi; mkdir -p /opt/longevidence-challenge/api'
rsync -rlpt -e "$SSH" api/server.py api/review_store.py api/review_views.py "$HOST:/opt/longevidence-challenge/api/"
rsync -rlpt -e "$SSH" deploy/longevidence-challenge-api.service "$HOST:/etc/systemd/system/longevidence-challenge-api.service"
rsync -rlpt -e "$SSH" deploy/add-api-route.py "$HOST:/root/add-challenge-api-route.py"
$SSH "$HOST" 'set -e
cat > /usr/local/bin/lc-review <<EOF
#!/bin/sh
# lc-review list | lc-review <id> accepted|rejected|unreviewed [note]
case "\$1" in list|review|credit) ;; *) set -- review "\$@" ;; esac
exec runuser -u lcapi -- env LC_DB=/var/lib/longevidence-challenge/challenge.db LC_TASKS_DIR=/opt/longevidence-challenge/challenge/tasks python3 /opt/longevidence-challenge/api/server.py "\$@"
EOF
chmod 755 /usr/local/bin/lc-review
systemctl daemon-reload
systemctl enable --now longevidence-challenge-api
systemctl restart longevidence-challenge-api
sleep 1
curl -sf http://127.0.0.1:8790/challenge/api/health
echo
python3 /root/add-challenge-api-route.py'
