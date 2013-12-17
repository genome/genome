-- Verify model.build_metric.index_build_id_metric_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'metric_build_id_metric_name_index';

ROLLBACK;
