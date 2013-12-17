-- Deploy web_nomenclature_field
-- requires: web_nomenclature

BEGIN;

CREATE TABLE IF NOT EXISTS web.nomenclature_field (
    id character varying(255) NOT NULL,
    name character varying(255) NOT NULL,
    type character varying(255) NOT NULL,
    nomenclature_id character varying(255) NOT NULL,
    CONSTRAINT nomenclature_field_pkey PRIMARY KEY (id),
    CONSTRAINT nomenclature_field_nomenclature_id_fkey FOREIGN KEY (nomenclature_id) REFERENCES web.nomenclature(id)
);

COMMIT;
