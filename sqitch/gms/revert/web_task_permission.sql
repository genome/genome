-- Revert web_task_permission

BEGIN;

REVOKE ALL ON TABLE web.task FROM "gms-user";

COMMIT;
