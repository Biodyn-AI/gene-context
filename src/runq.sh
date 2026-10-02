#!/bin/bash
# Serial job queue (the box has ~14 GB of real headroom and is shared): runs queue/todo/*.sh one at a time in name
# order, each after waiting for its "# mem N" GB of free+inactive memory. New jobs can be dropped into queue/todo
# while it runs. Exits when the queue is empty. Master log: queue/master.log; per-job logs: queue/logs/.
cd "$(dirname "$0")"; Q=queue; M=$Q/master.log
avail_gb() { vm_stat | awk '/page size of/ {ps=$8} /Pages free/ {f=$3} /Pages inactive/ {i=$3} /Pages speculative/ {s=$3} END {printf "%d", (f+i+s)*ps/1e9}'; }
wait_mem() { for k in $(seq 1 240); do a=$(avail_gb); [ "$a" -ge "$1" ] && return 0; echo "  [mem] ${a} GB < $1 GB; waiting" >> $M; sleep 30; done; return 0; }
echo "[runner start] $(date '+%F %T')" >> $M
while true; do
  next=$(ls $Q/todo/ 2>/dev/null | grep '\.sh$' | sort | head -1)
  [ -z "$next" ] && break
  mv "$Q/todo/$next" "$Q/running_$next"
  need=$(grep -m1 '^# mem' "$Q/running_$next" | awk '{print $3}'); need=${need:-4}
  wait_mem $need
  echo "[start] $next $(date '+%F %T') (avail $(avail_gb) GB)" >> $M
  bash "$Q/running_$next" > "$Q/logs/${next%.sh}.log" 2>&1; rc=$?
  echo "[end]   $next $(date '+%F %T') exit $rc" >> $M
  mv "$Q/running_$next" "$Q/done/$next"
done
echo "[queue empty] $(date '+%F %T')" >> $M
