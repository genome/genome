-- Deploy web_search_index_queue
-- requires: web_schema

BEGIN;

CREATE TABLE IF NOT EXISTS web.search_index_queue (
    id character varying(32) NOT NULL,
    subject_id character varying(256) NOT NULL,
    subject_class character varying(255) NOT NULL,
    "timestamp" timestamp(6) without time zone NOT NULL,
    priority smallint,
    CONSTRAINT siq_pk PRIMARY KEY (id)
);

COMMIT;
