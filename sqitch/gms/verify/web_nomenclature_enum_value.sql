-- Verify web_nomenclature_enum_value

BEGIN;

SELECT id, value, nomenclature_field_id
FROM web.nomenclature_enum_value
WHERE FALSE;

ROLLBACK;
