-- Deploy timeline_base
-- requires: timeline_schema

BEGIN;

CREATE TABLE IF NOT EXISTS timeline.base (
    id character varying NOT NULL,
    created_by character varying NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    name character varying NOT NULL,
    object_id character varying NOT NULL,
    object_class_name character varying NOT NULL,
    reason character varying NOT NULL,
    CONSTRAINT base_pkey PRIMARY KEY (id)
);

COMMIT;
