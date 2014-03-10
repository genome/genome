-- Deploy model_build_permission
-- requires: model_build

BEGIN;

REVOKE ALL ON TABLE model.build FROM PUBLIC;
REVOKE ALL ON TABLE model.build FROM genome;
GRANT ALL ON TABLE model.build TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.build TO "gms-user";

COMMIT;
