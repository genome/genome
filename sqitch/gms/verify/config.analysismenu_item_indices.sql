-- Verify config.analysismenu_item_indices

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'analysismenu_item_name_idx';
SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'analysismenu_item_path_idx';

ROLLBACK;
