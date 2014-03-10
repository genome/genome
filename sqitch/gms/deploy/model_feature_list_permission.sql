-- Deploy model_feature_list_permission
-- requires: model_feature_list

BEGIN;

REVOKE ALL ON TABLE model.feature_list FROM PUBLIC;
REVOKE ALL ON TABLE model.feature_list FROM genome;
GRANT ALL ON TABLE model.feature_list TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.feature_list TO "gms-user";

COMMIT;
