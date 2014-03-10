-- Deploy web.nomenclature.name
-- requires: web_nomenclature

BEGIN;

CREATE INDEX nomenclature_name_index on web.nomenclature using btree (name);

COMMIT;
