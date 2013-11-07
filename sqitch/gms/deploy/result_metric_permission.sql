-- Deploy result_metric_permission
-- requires: result_metric

BEGIN;

REVOKE ALL ON TABLE result.metric FROM PUBLIC;
REVOKE ALL ON TABLE result.metric FROM genome;
GRANT ALL ON TABLE result.metric TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE result.metric TO "gms-user";

COMMIT;
