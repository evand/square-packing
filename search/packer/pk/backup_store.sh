#!/usr/bin/env bash
# Daily store backup (10-15): online-safe sqlite snapshot (.backup works while chain / finisher write), zstd, keep the
# newest KEEP automatic copies.  Hand-made backups (any other name in the directory) are never touched.
# Installed as the systemd user timer pk-store-backup.timer (~/.config/systemd/user/); run by hand: pk/backup_store.sh
set -euo pipefail
STORE=${PK_STORE:-/home/evand/math/square-packing/public/search/packer/runs/store.sqlite}
DIR=${PK_BACKUP_DIR:-$HOME/store-backups}
KEEP=${PK_BACKUP_KEEP:-5}
mkdir -p "$DIR"
stamp=$(date +%Y-%m-%dT%H%M)
tmp="$DIR/.store-auto-$stamp.sqlite"
out="$DIR/store-auto-$stamp.sqlite.zst"
trap 'rm -f "$tmp"' EXIT
nice -n 19 sqlite3 "file:$STORE?mode=ro" ".timeout 600000" ".backup '$tmp'"
[ "$(sqlite3 "$tmp" 'PRAGMA quick_check')" = ok ] || { echo "backup: quick_check failed" >&2; exit 1; }
nice -n 19 zstd -q -T2 "$tmp" -o "$out"
rows=$(sqlite3 "$tmp" 'SELECT count(*) FROM packing')
echo "backup: $out ($(du -h "$out" | cut -f1), $rows packing rows)"
ls -1t "$DIR"/store-auto-*.sqlite.zst | tail -n +$((KEEP + 1)) | xargs -r rm -f --
