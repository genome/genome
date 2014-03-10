-- Deploy model.build_metric.build_id_metric_name
-- requires: model_build_metric

BEGIN;

CREATE INDEX metric_build_id_metric_name_index on model.build_metric using btree (build_id, metric_name);

COMMIT;
