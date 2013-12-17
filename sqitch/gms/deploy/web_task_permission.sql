-- Deploy web_task_permission
-- requires: web_task

BEGIN;

REVOKE ALL ON TABLE web.task FROM PUBLIC;
REVOKE ALL ON TABLE web.task FROM genome;
GRANT ALL ON TABLE web.task TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.task TO "gms-user";

COMMIT;
