-- Verify web_nomenclature

BEGIN;

SELECT id, name, default_value, accepts_any_field, empty_equivalent
FROM web.nomenclature
WHERE FALSE;

ROLLBACK;
