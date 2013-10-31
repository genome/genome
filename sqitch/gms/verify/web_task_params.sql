-- Verify web_task_params

BEGIN;

SELECT genome_task_id, params
FROM web.task_params
WHERE FALSE;

ROLLBACK;
