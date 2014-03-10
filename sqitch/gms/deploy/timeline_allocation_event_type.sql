-- Deploy timeline_allocation_event_type
-- requires: timeline_schema

BEGIN;

CREATE TABLE IF NOT EXISTS timeline.allocation_event_type (
    id character varying NOT NULL,
    CONSTRAINT allocation_event_type_pkey PRIMARY KEY (id)
);

COMMIT;
