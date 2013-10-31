-- Deploy subject_misc_attribute
-- requires: subject_schema

BEGIN;

CREATE TABLE IF NOT EXISTS subject.misc_attribute (
    entity_class_name character varying(255) NOT NULL,
    entity_id character varying(1000) NOT NULL,
    property_name character varying(255) NOT NULL,
    value character varying(4000),
    CONSTRAINT misc_attribute_pkey PRIMARY KEY (entity_id, entity_class_name, property_name)
);

COMMIT;
