-- Deploy web_nomenclature_enum_value
-- requires: web_nomenclature_field

BEGIN;

CREATE TABLE IF NOT EXISTS web.nomenclature_enum_value (
    id character varying(255) NOT NULL,
    value character varying(255) NOT NULL,
    nomenclature_field_id character varying(255) NOT NULL,
    CONSTRAINT nomenclature_enum_value_pkey PRIMARY KEY (id),
    CONSTRAINT nomenclature_enum_value_nomenclature_field_id_fkey
        FOREIGN KEY (nomenclature_field_id) REFERENCES web.nomenclature_field(id)
);

COMMIT;
