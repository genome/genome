-- Verify result.metric.index_metric_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_result_metric_metric_name';

ROLLBACK;
