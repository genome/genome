-- Deploy workflow_historian
-- requires: workflow_schema

BEGIN;

CREATE TABLE workflow.historian (
    net_key character varying(22) NOT NULL,
    operation_id numeric(19,0) NOT NULL,
    color   numeric(10,0) NOT NULL,
    workflow_instance_id character varying(255) NOT NULL,
    CONSTRAINT historian_pkey PRIMARY KEY (net_key, operation_id, color)
);

COMMIT;
