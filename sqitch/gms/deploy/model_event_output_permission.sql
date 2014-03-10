-- Deploy model_event_output_permission
-- requires: model_event_output

BEGIN;

REVOKE ALL ON TABLE model.event_output FROM PUBLIC;
REVOKE ALL ON TABLE model.event_output FROM genome;
GRANT ALL ON TABLE model.event_output TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.event_output TO "gms-user";

COMMIT;
