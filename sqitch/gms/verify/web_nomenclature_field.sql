-- Verify web_nomenclature_field

BEGIN;

SELECT id, name, type, nomenclature_id
FROM web.nomenclature_field
WHERE FALSE;

ROLLBACK;
