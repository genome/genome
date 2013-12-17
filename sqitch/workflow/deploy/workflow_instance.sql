-- Deploy workflow_instance
-- requires: workflow_plan
-- requires: workflow_instance_execution

BEGIN;

CREATE TABLE workflow.instance (
    workflow_instance_id character varying(255) NOT NULL,
    parent_instance_id character varying(255),
    peer_instance_id character varying(255),
    current_execution_id character varying(255),
    workflow_plan_id character varying(255) NOT NULL,
    name character varying(255) NOT NULL,
    input_stored bytea,
    output_stored bytea,
    parallel_index character varying(255),
    parent_execution_id character varying(255),
    intention character varying(15),
    CONSTRAINT instance_pkey PRIMARY KEY (workflow_instance_id),
    CONSTRAINT instance_peer_instance_id_fkey
        FOREIGN KEY (peer_instance_id) REFERENCES workflow.instance(workflow_instance_id),
    CONSTRAINT instance_parent_instance_id_fkey
        FOREIGN KEY (parent_instance_id) REFERENCES workflow.instance(workflow_instance_id),
    CONSTRAINT instance_workflow_plan_id_fkey
        FOREIGN KEY (workflow_plan_id) REFERENCES workflow.plan(workflow_plan_id),
    CONSTRAINT instance_current_execution_id_fkey
        FOREIGN KEY (current_execution_id) REFERENCES workflow.instance_execution(workflow_execution_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
