-- Deploy workflow_instance_execution
-- requires: workflow_schema

BEGIN;

CREATE TABLE workflow.instance_execution (
    workflow_execution_id character varying(255) NOT NULL,
    workflow_instance_id character varying(255) NOT NULL,
    status character varying(15),
    start_time timestamp without time zone,
    end_time timestamp without time zone,
    exit_code numeric(5,0),
    stdout character varying(512),
    stderr character varying(512),
    is_done numeric(2,0),
    is_running numeric(2,0),
    dispatch_id character varying(10),
    cpu_time numeric(15,0),
    max_memory numeric(10,0),
    max_swap numeric(10,0),
    max_processes numeric(4,0),
    max_threads numeric(4,0),
    user_name character varying(20),
    CONSTRAINT instance_execution_pkey PRIMARY KEY (workflow_execution_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
