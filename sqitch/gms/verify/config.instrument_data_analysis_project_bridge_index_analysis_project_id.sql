-- Verify config_analysis_menu_item_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instrument_data_analysis_project_bridge_analysis_project_id_idx';

ROLLBACK;
