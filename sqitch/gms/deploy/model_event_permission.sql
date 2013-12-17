-- Deploy model_event_permission
-- requires: model_event

BEGIN;

REVOKE ALL ON TABLE model.event FROM PUBLIC;
REVOKE ALL ON TABLE model.event FROM genome;
GRANT ALL ON TABLE model.event TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.event TO "gms-user";

COMMIT;
