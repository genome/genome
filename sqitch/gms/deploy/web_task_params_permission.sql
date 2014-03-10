-- Deploy web_task_params_permission
-- requires: web_task_params

BEGIN;

REVOKE ALL ON TABLE web.task_params FROM PUBLIC;
REVOKE ALL ON TABLE web.task_params FROM genome;
GRANT ALL ON TABLE web.task_params TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.task_params TO "gms-user";

COMMIT;
