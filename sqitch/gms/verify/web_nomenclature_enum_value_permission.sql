-- Verify web_nomenclature_enum_value_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.nomenclature_enum_value', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.nomenclature_enum_value', 'SELECT')::int;

ROLLBACK;
