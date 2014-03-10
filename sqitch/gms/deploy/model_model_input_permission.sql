-- Deploy model_model_input_permission
-- requires: model_model_input

BEGIN;

REVOKE ALL ON TABLE model.model_input FROM PUBLIC;
REVOKE ALL ON TABLE model.model_input FROM genome;
GRANT ALL ON TABLE model.model_input TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.model_input TO "gms-user";

COMMIT;
