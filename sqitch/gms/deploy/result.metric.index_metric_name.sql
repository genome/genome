-- Deploy result.metric.metric_name
-- requires: result_metric

BEGIN;

CREATE INDEX idx_result_metric_metric_name on result.metric using btree (metric_name);

COMMIT;
