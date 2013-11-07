-- Revert web_task_params_permission

BEGIN;

REVOKE ALL ON TABLE web.task_params FROM "gms-user";

COMMIT;
