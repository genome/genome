-- Verify config_analysis_menu_item_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'idx_allocation_path_like';

ROLLBACK;
