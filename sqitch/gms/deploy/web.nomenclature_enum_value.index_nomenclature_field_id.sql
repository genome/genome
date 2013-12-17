-- Deploy web.nomenclature_enum_value.nomenclature_field_id
-- requires: web_nomenclature_enum_value

BEGIN;

CREATE INDEX nomenclature_enum_field_index on web.nomenclature_enum_value using btree (nomenclature_field_id);

COMMIT;
