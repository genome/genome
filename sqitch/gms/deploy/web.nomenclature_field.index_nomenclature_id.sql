-- Deploy web.nomenclature_field.nomenclature_id
-- requires: web_nomenclature_field

BEGIN;

CREATE INDEX nomenclature_field_nomenclature_id_index on web.nomenclature_field using btree (nomenclature_id);

COMMIT;
