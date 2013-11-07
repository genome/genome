-- Verify web_nomenclature_field_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.nomenclature_field', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.nomenclature_field', 'SELECT')::int;

ROLLBACK;
