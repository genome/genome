-- Revert config.analysis_project.index_analysis_menu_item_id

BEGIN;

DROP INDEX config.c_ap_analysis_menu_item_id_index;

COMMIT;
