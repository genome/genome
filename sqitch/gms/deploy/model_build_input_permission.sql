-- Deploy model_build_input_permission
-- requires: model_build_input

BEGIN;

REVOKE ALL ON TABLE model.build_input FROM PUBLIC;
REVOKE ALL ON TABLE model.build_input FROM genome;
GRANT ALL ON TABLE model.build_input TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.build_input TO "gms-user";

COMMIT;
