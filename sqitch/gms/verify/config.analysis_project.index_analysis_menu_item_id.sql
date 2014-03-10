-- Verify config_analysis_menu_item_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'c_ap_analysis_menu_item_id_index';

ROLLBACK;
