-- Deploy workflow_execution_metric
-- requires: workflow_instance_execution

BEGIN;

CREATE TABLE workflow.execution_metric (
    workflow_execution_id character varying(255) NOT NULL,
    name character varying(255) NOT NULL,
    value character varying(1000),
    CONSTRAINT execution_metric_pkey PRIMARY KEY (workflow_execution_id, name),
    CONSTRAINT execution_metric_workflow_execution_id_fkey
        FOREIGN KEY (workflow_execution_id) REFERENCES workflow.instance_execution(workflow_execution_id)
);

COMMIT;
