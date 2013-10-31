-- Revert web_nomenclature

BEGIN;

DROP TABLE IF EXISTS web.nomenclature;

COMMIT;
