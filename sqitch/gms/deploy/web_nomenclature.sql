-- Deploy web_nomenclature
-- requires: web_schema

BEGIN;

CREATE TABLE IF NOT EXISTS web.nomenclature (
    id character varying(255) NOT NULL,
    name character varying(255) NOT NULL,
    default_value character varying(255),
    accepts_any_field boolean,
    empty_equivalent character varying(255),
    CONSTRAINT nomenclature_pkey PRIMARY KEY (id)
);

COMMIT;
