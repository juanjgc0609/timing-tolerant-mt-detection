#!/bin/bash
# vigila la Tarea 2L: relanza la cola si murio y no hay bandera de fin; imprime el progreso
cd /home/claude/work
if [ ! -f task2l_run_done.flag ] && [ "$(pgrep -fc 'python3 task2l_ru[n].py')" = "0" ]; then
  echo "relanzado $(date -u) (codigo sin cambios)" >> task2l_freeze_time.txt
  setsid nohup ./run_queue_2l.sh > /dev/null 2>&1 < /dev/null &
  sleep 5
fi
tail -1 task2l_run.log | cut -c1-120
ls task2l_run_done.flag 2>/dev/null; true
