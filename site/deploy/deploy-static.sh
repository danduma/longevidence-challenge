#!/bin/sh
# Publish dist/challenge to the existing Longevidence host.
# Backs up the live directory on the server first, then mirrors local files (removing obsolete ones).
# Usage: SSH_KEY=/path/to/key ./deploy/deploy-static.sh
set -eu
cd "$(dirname "$0")/.."
HOST=${HOST:?Set HOST to the deployment SSH destination}
SSH_KEY=${SSH_KEY:?Set SSH_KEY to your deployment key path}
SSH="ssh -o BatchMode=yes -o ConnectTimeout=10 -i $SSH_KEY"
STAMP=$(date -u +%Y%m%dT%H%M%SZ)

$SSH "$HOST" "cp -a /opt/longevidence-challenge/challenge /opt/longevidence-challenge/challenge.backup-$STAMP"
rsync -rlpt --delete -e "$SSH" dist/challenge/ "$HOST:/opt/longevidence-challenge/challenge/"
echo "Backup on server: /opt/longevidence-challenge/challenge.backup-$STAMP"

for path in /challenge/ /challenge/style.css /challenge/app.js /challenge/join.txt /challenge/tasks/lc-001.json /challenge/assets/archivo-expanded.woff2; do
  printf '%s %s\n' "$(curl -s -o /dev/null -w '%{http_code}' "https://longevidence.org$path")" "$path"
done
local_sha=$(shasum -a 256 dist/challenge/index.html | cut -d' ' -f1)
live_sha=$(curl -s https://longevidence.org/challenge/ | shasum -a 256 | cut -d' ' -f1)
[ "$local_sha" = "$live_sha" ] && echo "Live index.html matches local ($local_sha)" || { echo "MISMATCH: live $live_sha local $local_sha"; exit 1; }
