-- Deploy model_event_metric_permission
-- requires: model_event_metric

BEGIN;

REVOKE ALL ON TABLE model.event_metric FROM PUBLIC;
REVOKE ALL ON TABLE model.event_metric FROM genome;
GRANT ALL ON TABLE model.event_metric TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE model.event_metric TO "gms-user";

COMMIT;
