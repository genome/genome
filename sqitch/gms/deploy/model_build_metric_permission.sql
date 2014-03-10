-- Deploy model_build_metric_permission
-- requires: model_build_metric

BEGIN;

REVOKE ALL ON TABLE model.build_metric FROM PUBLIC;
REVOKE ALL ON TABLE model.build_metric FROM genome;
GRANT ALL ON TABLE model.build_metric TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.build_metric TO "gms-user";

COMMIT;
