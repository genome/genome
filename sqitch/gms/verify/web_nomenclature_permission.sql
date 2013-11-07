-- Verify web_nomenclature_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.nomenclature', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.nomenclature', 'SELECT')::int;

ROLLBACK;
