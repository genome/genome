-- Deploy model_event_input_permission
-- requires: model_event_input

BEGIN;

REVOKE ALL ON TABLE model.event_input FROM PUBLIC;
REVOKE ALL ON TABLE model.event_input FROM genome;
GRANT ALL ON TABLE model.event_input TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.event_input TO "gms-user";

COMMIT;
