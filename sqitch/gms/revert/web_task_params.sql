-- Revert web_task_params

BEGIN;

DROP TABLE IF EXISTS web.task_params;

COMMIT;
