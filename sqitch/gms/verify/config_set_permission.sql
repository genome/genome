-- Verify config_set_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'config.set', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'config.set', 'SELECT')::int;

ROLLBACK;
