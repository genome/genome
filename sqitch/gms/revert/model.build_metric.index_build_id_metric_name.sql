-- Revert model.build_metric.index_build_id_metric_name

BEGIN;

DROP INDEX model.metric_build_id_metric_name_index;

COMMIT;
