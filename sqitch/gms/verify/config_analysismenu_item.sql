-- Verify config_analysismenu_item

BEGIN;

SELECT *
FROM config.analysismenu_item
WHERE FALSE;

SELECT 1/has_table_privilege('genome', 'config.analysismenu_item', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'config.analysismenu_item', 'SELECT')::int;

ROLLBACK;
