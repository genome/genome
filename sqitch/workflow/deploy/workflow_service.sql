-- Deploy workflow_service
-- requires: workflow_schema

BEGIN;

CREATE TABLE workflow.service (
    hostname character varying(255) NOT NULL,
    port numeric(7,0) NOT NULL,
    start_time timestamp without time zone NOT NULL,
    process_id numeric(10,0) NOT NULL,
    username character varying(20) NOT NULL,
    CONSTRAINT service_pkey PRIMARY KEY (hostname, username, process_id, port, start_time)
);


COMMIT;
