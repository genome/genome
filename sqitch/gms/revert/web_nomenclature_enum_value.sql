-- Revert web_nomenclature_enum_value

BEGIN;

DROP TABLE IF EXISTS web.nomenclature_enum_value;

COMMIT;
