-- Revert result.metric.index_metric_name

BEGIN;

DROP INDEX result.idx_result_metric_metric_name;

COMMIT;
