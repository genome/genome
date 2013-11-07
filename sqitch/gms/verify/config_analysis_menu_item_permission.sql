-- Verify config_analysis_menu_item_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'config.analysis_menu_item', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'config.analysis_menu_item', 'SELECT')::int;

ROLLBACK;
